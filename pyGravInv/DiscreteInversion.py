import numpy as np
import math


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


def main():
    # load old inversion results
    depths, dens = np.loadtxt("../tests/data/output_inversion_L2.dens", unpack=True)
    measurement_depths = np.array((25., 66., 143.))
    for z_j in measurement_depths:
        Bj = row_of_B(z_j, depths)
        # with this method the data after free air gradient correction is returned, so to compare with the input
        # of the inversion, add the free air gradient back
        print(np.dot(Bj, dens) + 0.308*z_j)


if __name__ == '__main__':
    main()