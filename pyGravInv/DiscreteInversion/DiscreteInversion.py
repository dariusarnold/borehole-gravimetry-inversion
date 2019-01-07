import numpy as np

from pyGravInv.DiscreteInversion.components import matrix_A, misfit_squared
from pyGravInv.DiscreteInversion.helpers import get_inversion_depth_steps




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



