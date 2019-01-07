import numpy as np


def get_NL(fname):
    measurement_depths, _, measurement_errors = np.loadtxt(fname, unpack=True)
    inversion_depths = get_inversion_depth_steps(measurement_depths)
    # N: number of measurement points
    N = len(measurement_depths)
    # L : number of inversion depths
    L = len(inversion_depths)
    return N, L


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