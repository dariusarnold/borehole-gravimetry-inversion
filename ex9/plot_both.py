import numpy as np
import matplotlib.pyplot as plt


def main():
    """
    Show results of 9.2
    :return:
    """
    # load data from density model
    depth_dm, dens_dm = np.loadtxt("test_dm.dens", unpack=True)
    # load resulting density model from inversion
    depth_inv, dens_inv = np.loadtxt("test_inversion.dens", unpack=True)
    # convert to kg/m3
    dens_inv *= 1000

    fig, ax1 = plt.subplots()
    ax1.plot(depth_dm, dens_dm, color="red", label="Density model")
    ax2 = ax1.twinx()
    ax1.plot(depth_inv, dens_inv, color="blue", label="Inversion result")
    ax1.legend()
    ax1.set_ylabel("Density (kg/m³)")
    ax1.set_xlabel("Depth (m)")

    plt.show()


if __name__ == '__main__':
    main()