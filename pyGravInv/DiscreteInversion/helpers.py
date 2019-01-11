import numpy as np

from pyGravInv.DiscreteInversion.components import matrix_A


def get_N(fname):
    """
    Get number of measurement points N and number of inversion_depth_steps
    :param fname:
    :return:
    """
    measurement_depths, _, _ = np.loadtxt(fname, unpack=True)
    # N: number of measurement points
    N = len(measurement_depths)
    return N


def get_inversion_depth_steps(measurement_depths, num_steps=10000):
    """
    Get the discretization depths of the inversion result
    :param measurement_depths: np array containing measurement depths
    :param num_steps: Number of steps used for discretozation
    :return:
    """
    return np.linspace(0, measurement_depths[-1], num_steps)


def correct_free_air_gradient(measurement_depths, measurement_data):
    """
    Correct given measurement data for the free air gradient
    :return: corrected measurement data
    """
    f = 0.308 # mGal / m
    return measurement_data - f * measurement_depths


def read_data(fname):
    """
    Read measurement data with errors from fname and correct the data for the free air gradient
    :param fname:
    :return: measurement depth, measurement data (free air gradient corrected), measurement errors
    """
    measurement_depths, measurement_data, measurement_errors = np.loadtxt(fname, unpack=True)
    measurement_data = correct_free_air_gradient(measurement_depths, measurement_data)
    return measurement_depths, measurement_data, measurement_errors


def matrix_A_from_file(fname, num_steps=10000):
    # read/calculate required values
    measurement_depths, _, measurement_errors = np.loadtxt(fname, unpack=True)
    inversion_depths = get_inversion_depth_steps(measurement_depths, num_steps=num_steps)
    # create Matrix A
    A = matrix_A(measurement_depths, measurement_errors, inversion_depths)
    return A


def do_SVD_from_file(fname, full_matrices=False, num_steps=10000):
    """

    :param fname:
    :return:
    """
    # create Matrix A
    A = matrix_A_from_file(fname, num_steps)
    # do svd
    V, Lambda, U_transposed = np.linalg.svd(A, full_matrices=full_matrices)
    return V, Lambda, U_transposed, num_steps