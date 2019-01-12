import matplotlib.pyplot as plt

from pyGravInv.DiscreteInversion.DiscreteInversion import resolution_matrix
from pyGravInv.DiscreteInversion.helpers import get_inversion_depth_steps, read_data


def ex_1(fname):
    NUMBER_OF_DEPTH_STEPS = 200
    measurement_depth, *_ = read_data(fname)
    inversion_depths = get_inversion_depth_steps(measurement_depth, num_steps=NUMBER_OF_DEPTH_STEPS)
    R = resolution_matrix(fname, depth_num_steps=NUMBER_OF_DEPTH_STEPS)
    extent = (inversion_depths[0], inversion_depths[-1], inversion_depths[-1], inversion_depths[0])
    plt.imshow(R, interpolation="none", extent=extent)
    cbar = plt.colorbar()
    cbar.set_label("Density (g/cm³)")
    plt.xlabel("Tiefe der dünnen Schicht (m)")
    plt.ylabel("Tiefe (m)")
    plt.title("Resolutionsmatrix")
    plt.show()


if __name__ == '__main__':
    ex_1("data/grav15.dat")
