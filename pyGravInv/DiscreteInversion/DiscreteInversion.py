import numpy as np

from pyGravInv.DiscreteInversion.components import new_data, calculate_coeffs, optimal_nu_bysection, \
    matrix_S_inverted_helper, calculate_model
from pyGravInv.DiscreteInversion.helpers import get_inversion_depth_steps, read_data, do_SVD_from_file
from pyGravInv.ModelMaker.DensityModel import DiscretizedDensityModel


class InversionModel(DiscretizedDensityModel):
    def __init__(self, depths, densities, nu, misfit):
        """
        Holds depth, density values of the result of an inversion
        :param depths: ndarray of depths values increasing from 0 to the last measurement depths
        :param densities: ndarray of density values associated with the depth values
        :param nu: Lagrange parameter nu
        :param misfit: Misfit of the density model.
        """
        super().__init__(depths, densities)
        self.nu = nu
        self.misfit = misfit

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
    """
    Do an inversion using the borehole gravity measurements with errors given in fname
    :param fname: Filepath to measurement file which contains three values per line
    depth (m) measurement data (mGal, not corrected for free air gradient), measurement error (mgal)
    Every line represents one measurement point. Measurement points in the file are ordered with increasing depth.
    :param nu: Lagrange parameter
    :param depth_num_steps: Number of steps for the resulting discrete density Model
    :return: InversionModel which holds the depth, density values of the resulting model
    """
    measurement_depths, measurement_data, measurement_errors = read_data(fname)
    inversion_depths = get_inversion_depth_steps(measurement_depths, num_steps=depth_num_steps)
    V, Lambda, U_transposed, _ = do_SVD_from_file(fname, num_steps=depth_num_steps)
    # calculate new data
    d_double_prime = new_data(measurement_data, np.transpose(V),  measurement_errors)
    # calculate desired misfit using bisection
    desired_misfit = len(measurement_depths)
    optimal_nu, misfit = optimal_nu_bysection(desired_misfit, d_double_prime, Lambda)
    # calculate new coefficients
    alpha_double_prime = calculate_coeffs(d_double_prime, optimal_nu, Lambda)
    S_inverted = matrix_S_inverted_helper(inversion_depths)
    dens_model = calculate_model(S_inverted, np.transpose(U_transposed), new_coefficients=alpha_double_prime)
    return InversionModel(inversion_depths, dens_model, optimal_nu, misfit)


def resolution_matrix(fname, nu=None, depth_num_steps=10000):
    """
    Calculate the resolution matrix R. Every column of the matrix is the inversion result for a test model containing
    one spike at a certain depth.
    :param fname: Filepath to measurement file which contains three values per line
    depth (m) measurement data (mGal, not corrected for free air gradient), measurement error (mgal)
    :param nu: Lagrange parameter
    :param depth_num_steps: Number of steps for the resulting discrete density Model
    :return: matrix R
    """
    measurement_depths, measurement_data, measurement_errors = read_data(fname)
    inversion_depths = get_inversion_depth_steps(measurement_depths, num_steps=depth_num_steps)
    S_inverted = matrix_S_inverted_helper(inversion_depths)
    V, Lambda, U_transposed, _ = do_SVD_from_file(fname, num_steps=depth_num_steps)
    # calculate new data
    d_double_prime = new_data(measurement_data, np.transpose(V),  measurement_errors)
    # calculate desired misfit using bisection
    desired_misfit = len(measurement_depths)
    optimal_nu, misfit = optimal_nu_bysection(desired_misfit, d_double_prime, Lambda)
    D_nu = matrix_D_nu(optimal_nu, Lambda)
    D_nu_inverted = np.linalg.inv(D_nu)
    R = matrix_R(S_inverted, U_transposed.T, D_nu_inverted, Lambda)
    return R
