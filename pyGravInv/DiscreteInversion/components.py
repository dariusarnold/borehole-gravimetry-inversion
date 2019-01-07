"""
Contains functions calculating all the required variables
"""
import numpy as np


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


def calculate_coeffs(new_data, nu, Lambda):
    return new_data / ( (1/nu * 1/Lambda) + Lambda )


def calculate_model(S_inverted, U, new_coefficients):
    return S_inverted @ U @ new_coefficients


def do_SVD(A, full_matrices=False):
    V, Lambda, U_transposed = np.linalg.svd(A, full_matrices=full_matrices)
    return V, Lambda, U_transposed


def optimal_nu_bysection(target_misfit, new_data, Lambda, accuracy=0.01):
    """
    Calculate the optimal Lagrange parameter using a bisection method.
    The optimum nu maximizes the misfit.
    :param target_misfit: The misfit to reach
    :return: Nu required to reach target misfit and the actual misfit given the accuracy of nu
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
    return nu_mid, chi_squared