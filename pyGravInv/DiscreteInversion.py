import numpy as np
import math
import pyGravInv


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


def ex11_1():
    # invert the data
    (inversion_depths, inversion_dens), _ = pyGravInv.inversion_error("data/grav15.dat", pyGravInv.ErrorNorm.L2ErrorNorm)
    measurement_depths, _, measurement_errors = np.loadtxt("data/grav15.dat", unpack=True)

    # create Matrix A
    A = matrix_A(measurement_depths, measurement_errors, inversion_depths)
    # do svd
    V, Lambda, U_transposed = np.linalg.svd(A, full_matrices=False)
    # N: number of measurement points
    N = len(measurement_depths)
    # L : number of inversion depths
    L = len(inversion_depths)

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




if __name__ == '__main__':
    ex11_1()