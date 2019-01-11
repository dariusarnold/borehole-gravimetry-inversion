import matplotlib.pyplot as plt
import numpy as np

from pyGravInv.DiscreteInversion.DiscreteInversion import invert_errors
from pyGravInv.DiscreteInversion.components import new_data, misfit_squared, matrix_S_inverted_helper
from pyGravInv.DiscreteInversion.helpers import get_N, read_data, get_inversion_depth_steps, do_SVD_from_file


def print_check_result(name, check):
    print_helper = lambda condition: "OK" if condition else "FAIL"
    s = f"{name:<14} {print_helper(check):.>6}"
    print(s)


def ex11_1(fname):
    """
    Do SVD (A = V S U.T) and check the following conditions:
    V.T V = V V.T = I_N
    U.T U = I_N
    U U.T != I_L
    for the NxL matrix A
    """
    # do just 1000 discretiation steps to get quicker results
    V, Lambda, U_transposed, L = do_SVD_from_file(fname, num_steps=1000)
    N = get_N(fname)

    # check if the conditions are fulfilled
    V_transposed = np.transpose(V)
    VTV = V_transposed @ V
    VVT = V @ V_transposed
    I_N = np.identity(N)
    I_L = np.identity(L)
    print_check_result(name="V^T*V == V*V^T", check=np.allclose(VTV, VVT))
    print_check_result(name="V^T*V == I_N", check=np.allclose(VTV, I_N))
    # next condition
    U = np.transpose(U_transposed)
    UUT = U @ U_transposed
    UTU = U_transposed @ U
    print_check_result(name="U^T*U == I_N", check=np.allclose(UTU, I_N))
    print_check_result(name="U*U^T != I_L", check=not np.allclose(UUT, I_L))


def ex11_2(fname):
    """
    Calculate new data d'' and misfit chi² and plot chi²(nu)
    """
    measurement_depths, measurement_data, measurement_errors = read_data(fname)
    V, Lambda, U_transposed, _ = do_SVD_from_file(fname)
    # calculate new data
    d_double_prime = new_data(measurement_data, np.transpose(V),  measurement_errors)
    # calculate misfit for multiple values of nu
    range_of_nus = np.logspace(-2, 2, 1000)
    range_of_chis = np.array([misfit_squared(d_double_prime, nu, Lambda) for nu in range_of_nus])
    plt.plot(range_of_nus, range_of_chis, ".")
    plt.title(r"$\nu/\chi^2$")
    plt.xlabel(r"Lagrangeparameter $\nu$")
    plt.ylabel(r"Misfit $\chi^2$")
    plt.xscale("log")
    plt.yscale("log")
    plt.show()


def ex11_3(fname):
    """
    Do inversion and plot result
    """
    depths, densities = invert_errors(fname)
    plt.plot(depths, densities)
    plt.title("Density Model")
    plt.xlabel(r"Depth $z$ (m)")
    plt.ylabel(r"Density $\rho$ (g/cm³)")
    plt.show()


def ex11_4(fname):
    """
    Plot eigenvectors of the null space and the orthogonal representers
    """
    measurement_depths, measurement_data, measurement_errors = read_data(fname)
    inversion_depths = get_inversion_depth_steps(measurement_depths, num_steps=200)
    V, Lambda, U_transposed, _ = do_SVD_from_file(fname, full_matrices=True, num_steps=200)
    S_inverted = matrix_S_inverted_helper(inversion_depths)
    U = np.transpose(U_transposed)
    representants = S_inverted @ U
    # representants now holds the representants in its columns, but it is set up so that indexing it accesses its rows
    # transpose it to access the representants by index
    representants = np.transpose(representants)

    # create enough subplots to plot all representes
    f, ax_arr = plt.subplots(len(measurement_data), sharex=True)
    plt.suptitle("Representants")
    for i, repr in enumerate(representants[0:len(measurement_data)]):
        # show singular value belonging to representer in legend
        ax_arr[i].plot(inversion_depths, repr, label=f"$\lambda_{{{i+1}}}$: {Lambda[i]:4.1f}")
        ax_arr[i].legend()
    plt.xlabel(r"Depth $z$ in m")
    plt.figure()
    for nullspace_vector in representants[len(measurement_data):]:
        plt.plot(inversion_depths, nullspace_vector)
    plt.xlabel(r"Depth $z$ in m")
    plt.title("Eigenvectors null space")
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
