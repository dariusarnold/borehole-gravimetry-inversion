import numpy as np
import math
import pyGravInv
import matplotlib.pyplot as plt


def find_closest_index(array, value):
    """
    Return the index of the elements that are closest to value in array
    :param array:
    :param value:
    :return:
    """
    return np.argmin(np.abs(array-value))


def find_closest_le(array, value):
    """
    Return the index of the element in array that is closest to value, but <= value
    :param array:
    :param value:
    :return:
    """
    return find_closest_index(array, math.floor(value))


def row_of_B(z_j, z_nu):
    """
    Calculate one row of the matrix B
    :param z_j: Scalar value up to which
    :type z_j: float
    :param z_nu: vector with grid points
    :type z_nu: np.ndarray
    :return:
    """

    def delta_z(index):
        """
        Calculate delta z for a given index as a forward difference
        :param index: index nu
        """
        return z_nu[index+1] - z_nu[index]

    alpha = find_closest_le(z_nu, z_j)
    phi = (z_j - z_nu[alpha]) / (z_nu[alpha+1] - z_nu[alpha])
    gamma = 0.08382 # mGal*cm³/(m*g)
    B = np.zeros(len(z_nu))
    for nu in range(0, len(z_nu)):
        # all indices in the conditions - 1 since in math the indices start at 1, but in python they start at 0
        if nu == 0:
            B[nu] = delta_z(nu)
        elif 0 < nu < alpha-1:
            B[nu] = delta_z(nu-1) + delta_z(nu)
        elif nu == alpha-1:
            B[nu] = delta_z(nu-1) + phi * delta_z(nu) * (2 - phi)
        elif nu == alpha:
            B[nu] = phi**2 * delta_z(nu-1)
        else:
            # B is default initialized as zero, so we dont have to set it here
            continue
    B *= -gamma/2
    return B


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
    Create full matrix B for a gravimetry inversion
    :param measurement_depths: numpy array containing the depths where values were measured
    :param inversion_depths: numpy array containing the depths where the inversion result was discretized
    :return:
    """
    # Create empty array. Since we want to append rows, we have to create it with the number of columns the rows are
    # going to have
    B = np.empty((0, len(inversion_depths)))
    for z in measurement_depths:
        B_row = row_of_B(z, inversion_depths)
        # append row to matrix B
        B = np.vstack((B, B_row))
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
    A = np.dot(Sigma_inverted, B)
    A = np.dot(A,  S_inverted)
    return A


def print_helper(condition):
    return "OK" if condition else "FAIL"


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


def ex11_1(fname):
    V, Lambda, U_transposed = do_SVD(fname)
    N, L = get_NL(fname)

    # check if the conditions are fullfilled
    V_transposed = np.transpose(V)
    VTV = np.dot(V_transposed, V)
    VVT = np.dot(V, V_transposed)
    I_N = np.identity(N)
    I_L = np.identity(L)
    print("V^T*V == V*V^T:", print_helper(np.allclose(VTV, VVT)))
    print("V^T*V == I_N:", print_helper(np.allclose(VTV, I_N)))
    #next condition
    U = np.transpose(U_transposed)
    UUT = np.dot(U, U_transposed)
    UTU = np.dot(U_transposed, U)
    print("U^TU == I_N:", print_helper(np.allclose(UTU, I_N)))
    print("UU^T != I_L:", print_helper(not np.allclose(UUT, I_L)))


def get_inversion_depth_steps(measurement_depths, num_steps=10000):
    """
    Get the discretization depths of the inversion result
    :param measurement_depths: np array containing measurement depths
    :param num_steps: Number of steps used for discretozation
    :return:
    """
    num_steps += 2
    return np.linspace(0, measurement_depths[-1], num_steps)


def new_data(old_data, V_transposed, measurement_errors):
    return (1/measurement_errors) * np.dot(V_transposed, old_data)


def misfit_squared(new_data, nu, Lambda):
    """
    Calculate new misfit
    :param new_data: d''
    :param nu: Lagrange parameter
    :param Lambda: Vector containing the Eigenvalues of the SVD
    :return: Chi²
    """
    zwischen = (new_data / (1 + nu*Lambda))**2
    return sum(zwischen)


def ex11_2(fname):
    measurement_depths, measurement_dens, measurement_errors = np.loadtxt(fname, unpack=True)
    V, Lambda, U_transposed = do_SVD(fname)
    # calculate new data
    d_double_prime = new_data(measurement_dens, np.transpose(V),  measurement_errors)
    # calculate misfit for multiple values of nu
    range_of_nus = np.logspace(-2, 2, 1000)
    range_of_chis = np.array([misfit_squared(d_double_prime, nu, Lambda) for nu in range_of_nus])
    plt.plot(range_of_nus, range_of_chis, ".")
    plt.show()


def calculate_coeffs(new_data, nu, Lambda):
    return new_data / ( (1/nu + 1/Lambda) + Lambda )


def calculate_model(S_inverted, U, new_coefficients):
    return np.dot(np.dot(S_inverted, U), new_coefficients)


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
    print(nu_mid)
    return nu_mid


def ex11_3(fname):
    measurement_depths, measurement_dens, measurement_errors = np.loadtxt(fname, unpack=True)
    inversion_depths = get_inversion_depth_steps(measurement_depths)
    V, Lambda, U_transposed = do_SVD(fname)
    desired_misfit = len(measurement_depths)
    # calculate new data
    d_double_prime = new_data(measurement_dens, np.transpose(V),  measurement_errors)
    # calculate desired misfit using bisection
    optimal_nu = optimal_nu_bysection(desired_misfit, d_double_prime, Lambda)
    # calculate new coefficients
    alpha_double_prime = calculate_coeffs(d_double_prime, optimal_nu, Lambda)
    S_shape = (len(inversion_depths), len(inversion_depths))
    dx = inversion_depths[1] - inversion_depths[0]
    S_inverted = matrix_S_inverted(gamma=1E-3, shape=S_shape, delta_x=dx)
    dens_model = calculate_model(S_inverted, np.transpose(U_transposed), new_coefficients=alpha_double_prime)
    plt.plot(dens_model)
    plt.show()



def main():
    #ex11_1("data/grav15.dat")
    ex11_2("data/grav15.dat")
    ex11_3("data/grav15.dat")


if __name__ == '__main__':
    main()