import numpy as np

from pyGravInv.DiscreteInversion.components import matrix_A, misfit_squared, new_data, calculate_coeffs, matrix_S_inverted, calculate_model
from pyGravInv.DiscreteInversion.helpers import get_inversion_depth_steps, read_data

class InversionModel:
    def __init__(self, depths, densities):
        """
        Holds depth, density values of the result of an inversion
        :param depths: ndarray of depths values increasing from 0 to the last measurement depths
        :param densities: ndarray of density values associated with the depth values
        """
        self.depths = depths
        self.densities = densities

    def __iter__(self):
        """
        Make the class work with tuple unpacking, eg:
        >>> i = InversionModel(1, 2)
        >>> depth, dens = i
        >>> depth
        1
        >>> dens
        2
        """
        return iter((self.depths, self.densities))


def invert_errors(fname, nu=None, depth_num_steps=10000):
    measurement_depths, measurement_data, measurement_errors = read_data(fname)
    inversion_depths = get_inversion_depth_steps(measurement_depths, num_steps=depth_num_steps)
    V, Lambda, U_transposed = do_SVD(fname, num_steps=depth_num_steps)
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
    return InversionModel(inversion_depths, dens_model)


#TODO write version that directly takes matrix A
def do_SVD(fname, full_matrices=False, num_steps=10000):
    """

    :param fname:
    :return:
    """
    # read/calculate required values
    measurement_depths, _, measurement_errors = np.loadtxt(fname, unpack=True)
    inversion_depths = get_inversion_depth_steps(measurement_depths, num_steps=num_steps)

    # create Matrix A
    A = matrix_A(measurement_depths, measurement_errors, inversion_depths)
    # do svd
    V, Lambda, U_transposed = np.linalg.svd(A, full_matrices=full_matrices)
    return V, Lambda, U_transposed


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
